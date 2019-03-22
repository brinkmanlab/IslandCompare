import * as Common from "./_common";
import { Model as VuexModel } from '@vuex-orm/core';
import { Model as History } from './histories';

class HistoryDatasetAssociation extends VuexModel {
    static entity = 'HistoryDatasetAssociation';
    static primaryKey = 'id';

    static fields() {
        return {
            accessible: this.boolean(false),
            type_id: this.string(null).nullable(),
            file_name: this.string(null).nullable(),
            resubmitted: this.boolean(false),
            create_time: this.string((new Date).toISOString()),
            creating_job: this.string(null).nullable(),
            dataset_id: this.string(null).nullable(),
            file_size: this.number(0),
            file_ext: this.string(null).nullable(),
            id: this.string(null).nullable(),
            misc_info: this.string(''),
            hda_ldda: this.string('hda'),
            download_url: this.string(null).nullable(),
            state: this.string('ok'),
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
        }
    }

    static methodConf = {
        http: {
            url: ':contents_url/datasets' //TODO
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '',
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

class HistoryDatasetCollectionAssociation extends VuexModel {
    static entity = 'HistoryDatasetCollectionAssociation';
    static primaryKey = 'id';

    static fields() {
        return {
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

    static methodConf = {
        http: {
            url: ':contents_url/dataset_collections'
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '',
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