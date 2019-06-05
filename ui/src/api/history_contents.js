import * as Common from "./_common";
import { History } from './histories';
import axios from "axios";
import { galaxy_path } from "@/app.config";

class HistoryDatasetAssociation extends Common.Model {
    static entity = 'HistoryDatasetAssociation';
    static primaryKey = 'id';
    static end_states = ['ok', 'error'];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            accessible: this.boolean(false),
            type_id: this.string(null).nullable(),
            file_name: this.string(null).nullable(),
            resubmitted: this.boolean(false),
            create_time: this.string((new Date).toISOString()),
            creating_job: this.string(null).nullable(),
            dataset_id: this.string(null).nullable(),
            file_size: this.number(0),
            file_ext: this.string(null).nullable(),
            misc_info: this.string(''),
            hda_ldda: this.string('hda'),
            download_url: this.string(null).nullable(),
            state: this.string('pending'),
            display_types: this.attr([]),
            display_apps: this.attr([]),
            metadata_dbkey: this.string('?'),
            type: this.string(null).nullable(),
            misc_blurb: this.string(''),
            peek: this.string(''),
            update_time: this.string((new Date).toISOString()),
            data_type: this.string(null).nullable(),
            tags: this.attr([]),
            deleted: this.boolean(false),
            history_id: this.string(null).nullable(),
            meta_files: this.attr([]),
            genome_build: this.string('?'),
            metadata_sequences: this.number(0),
            hid: this.number(0),
            model_class: this.string('HistoryDatasetAssociation'),
            metadata_data_lines: this.number(0),
            annotation: this.string(null).nullable(),
            permissions: this.attr({
                access: [],
                manage: [],
            }),
            history_content_type: this.string('dataset'),
            name: this.string('New dataset'),
            extension: this.string(''),
            visible: this.boolean(true),
            url: this.string(null).nullable(),
            uuid: this.string(null).nullable(),
            visualizations: this.attr([]),
            rerunnable: this.boolean(false),
            purged: this.boolean(false),
            api_type: this.string(null).nullable(),

            //ORM only
            history: this.belongsTo(History, 'history_id'),
            //TODO dataset: this.hasOne(Dataset, 'dataset_id'),
            file: this.attr(null),
            upload_progress: this.number(100),
        }
    }

    get_base_url() {
        let history = this.history;
        if (!history) history = History.find(this.history_id);
        return history.contents_url;
    }

    //TODO move to Dataset, this was placed here because the api returns an hda object
    static async $upload(file, history_id, http) { //eslint-disable-line
        let tmp_id = file.name+Math.floor(Math.random()*10**16).toString();
        await HistoryDatasetAssociation.insert({data: {id: tmp_id, file: file, name: file.name, hid: 0, history_id: history_id, upload_progress: 0, state: "uploading"}});
        //let payload = {
        //    tool_id: 'upload1',
        //    history_id: this.history_id,
        //    inputs: {
        //        'files_0|NAME'  : this.file_name,
        //        'files_0|type'  : 'upload_dataset',
        //        'dbkey'         : '?',
        //        'file_type'     : 'auto',
        //        'ajax_upload'   : 'true',
        //    },
        //    'files_0|file_data': file,
        //};
        //let data = Object.entries(payload).reduce((form, [key, val])=>{ form.append(key, (val instanceof Object) ? JSON.stringify(val) : val); return form; }, new FormData());
        let formData = new FormData();
        formData.append('history_id', history_id);
        let inputs = {
            'files_0|NAME'  : file.name,
            'files_0|type'  : 'upload_dataset',
            'dbkey'         : '?',
            'file_type'     : 'auto',
            'ajax_upload'   : 'true',
        };
        formData.append('inputs', JSON.stringify(inputs));
        formData.append('tool_id', 'upload1');
        formData.append('files_0|file_data', file);
        let response = await axios.post('/api/tools/', formData, {
            baseURL: galaxy_path,
            params: {
                //key: 'admin',
                uuid: this.store().state.user_uuid, //TODO why store() and not $store
            },
            headers: {
                'Content-Type': 'application/json',
            },
            onUploadProgress: progressEvent => {
                //console.log(progressEvent.total); // eslint-disable-line no-console
                let upload_progress = (progressEvent.loaded * 100 / progressEvent.total) - 1; //-1 To keep progress indicator active until temporary hda deleted
                HistoryDatasetAssociation.update({
                    where: tmp_id,
                    data: { upload_progress: upload_progress }
                });
            },
            ...http,
        });
        //let response = await Tool.$create({
        //    /*http: {
        //        headers: {'Content-Type': 'multipart/form-data'},
        //        //...http,
        //    },*/
        //    data: formData,
        //    //TODO transformResponse to extract HDA
        //});
        if (response.status !== 200) {
            HistoryDatasetAssociation.update({
                where: tmp_id,
                data: { name: "Upload failed" } //TODO append file name
            });
        } else {
            HistoryDatasetAssociation.delete(tmp_id);
            return await HistoryDatasetAssociation.insert({data: response.data.outputs[0]});
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: ':url/datasets',
            /*transformRequest: [function transformRequest(data, headers) {
                console.log(data, headers); //eslint-disable-line
                //if (data instanceof FormData) {
                //    headers['Content-Type'] = 'multipart/form-data';
                //}
                return data;
            }],*/
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '?view=detailed',
                    method: 'get',
                },
            },
            $get: {
                name: 'get',
                http: {
                    url: '/:id',
                    method: 'get',
                },
            },
            $create: {
                name: 'create',
                http: {
                    url: '', //TODO
                    method: 'post',
                },
            },
            $update: {
                name: 'update',
                http: {
                    url: '/:id', //TODO
                    method: 'put',
                },
            },
            $delete: {
                name: 'delete',
                http: {
                    url: '/:id', //TODO
                    method: 'delete',
                },
            },
        }
    }
}

class HistoryDatasetCollectionAssociation extends Common.Model {
    static entity = 'HistoryDatasetCollectionAssociation';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            //Request data
            id: this.string(null).nullable(),
            history_content_type: this.string('dataset_collection'),
            populated_state_message: this.string(null).nullable(),
            name: this.string('New collection'),
            populated: this.boolean(false),
            deleted: this.boolean(false),
            type: this.string('collection'),
            history_id: this.string(null).nullable(),
            tags: this.attr([]),
            visible: this.boolean(true),
            job_source_id: this.string(null).nullable(), //TODO add relation
            job_source_type: this.string(null).nullable(),
            collection_type: this.string('list'),
            url: this.string(null).nullable(),
            model_class: this.string('HistoryDatasetCollectionAssociation'),
            hid: this.number(0),
            element_count: this.number(0),
            populated_state: this.string('ok'),
            elements: this.attr([]), //TODO add relation

            //ORM only
            history: this.belongsTo(History, 'history_id'),
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: ':url/dataset_collections'
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '?view=detailed',
                    method: 'get',
                },
            },
            $get: {
                name: 'get',
                http: {
                    url: '/:id',
                    method: 'get',
                },
            },
            $create: {
                name: 'create',
                http: {
                    url: '', //TODO
                    method: 'post',
                },
            },
            $update: {
                name: 'update',
                http: {
                    url: '/:id', //TODO
                    method: 'put',
                },
            },
            $delete: {
                name: 'delete',
                http: {
                    url: '/:id', //TODO
                    method: 'delete',
                },
            },
        }
    }
}

const Module = {
    ...Common.Module,
    state: Common.State({

    }),
    mutations: {
        ...Common.Mutations,
    },
    actions: {
        ...Common.Actions,
    },
    getters: {
        ...Common.Getters,
    },
};

function register(database) {
    database.register(HistoryDatasetAssociation, Module);
    database.register(HistoryDatasetCollectionAssociation, Module);
}

export {
    HistoryDatasetAssociation,
    HistoryDatasetCollectionAssociation,
    Module,
    register,
};