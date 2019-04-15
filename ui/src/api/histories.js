import * as Common from "./_common";

import { HistoryDatasetAssociation, HistoryDatasetCollectionAssociation } from "./history_contents"; //eslint-disable-line


class History extends Common.Model {
    static entity = 'History';
    static primaryKey = 'id';
    static end_states = ['ok','error'];

    static fields() {
        return {
            ...super.fields(),
            //Request data
            id: this.string(null),
            importable: this.boolean(false),
            create_time: this.string((new Date).toISOString()),
            contents_url: this.string(null).nullable(),
            size: this.number(0),
            user_id: this.string(null).nullable(),
            username_and_slug: this.string(null).nullable(),
            annotation: this.string(null).nullable(),
            state_details: this.attr({
                paused: 0,
                ok: 0,
                failed_metadata: 0,
                upload: 0,
                discarded: 0,
                running: 0,
                setting_metadata: 0,
                error: 0,
                new: 0,
                queued: 0,
                empty: 0,
            }),
            state: this.string('ok'),
            empty: this.boolean(true),
            update_time: this.string((new Date).toISOString()),
            tags: this.attr([]),
            deleted: this.boolean(false),
            genome_build: this.string(null).nullable(),
            slug: this.string(null).nullable(),
            name: this.string("Unnamed History"),
            url: this.string(null).nullable(),
            state_ids: this.attr({
                paused: [],
                ok: [],
                failed_metadata: [],
                upload: [],
                discarded: [],
                running: [],
                setting_metadata: [],
                error: [],
                new: [],
                queued: [],
                empty: [],
            }),
            published: this.boolean(false),
            model_class: this.string("History"),
            purged: this.boolean(false),

            //ORM only
            datasets: this.hasMany(HistoryDatasetAssociation, 'history_id'),
            collections: this.hasMany(HistoryDatasetCollectionAssociation, 'history_id'),
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/histories'
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
                    url: '?view=detailed',
                    method: 'post',
                },
            },
            $update: {
                name: 'update',
                http: {
                    url: '/:id',
                    method: 'put',
                },
            },
            $delete: {
                name: 'delete',
                http: {
                    url: '/:id',
                    method: 'delete',
                },
            },
        }
    }
}

/*class HistoryExport extends Common.Model { //TODO
    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/histories'
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
                    url: '?view=detailed',
                    method: 'post',
                },
            },
            $update: {
                name: 'update',
                http: {
                    url: '/:id',
                    method: 'put',
                },
            },
            $delete: {
                name: 'delete',
                http: {
                    url: '/:id',
                    method: 'delete',
                },
            },
        }
    }
}*/

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
    database.register(History, Module);
}

export {
    History,
    Module,
    register,
};