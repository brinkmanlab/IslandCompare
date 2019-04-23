import * as Common from "./_common";
//import axios from "axios";
//import {HistoryDatasetAssociation, HistoryDatasetCollectionAssociation } from "./history_contents";
import { History } from "./histories";


class WorkflowInvocationStep extends Common.Model {
    static entity = 'WorkflowInvocationStep';
    static primaryKey = 'id';
    static end_states = ['scheduled','error'];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            workflow_step_uuid: this.string(null).nullable(),
            update_time: this.string(null).nullable(),
            jobs: this.attr([]),
            job_id: this.string(null).nullable(),
            outputs: this.attr(), //this.hasMany(HistoryDatasetAssociation, 'id'), //TODO
            order_index: this.number(0),
            output_collections: this.attr(), //this.hasMany(HistoryDatasetCollectionAssociation, 'id'), //TODO
            workflow_step_label: this.string(''),
            state: this.string(null).nullable(),
            action: this.string(null).nullable(),
            model_class: this.string("WorkflowInvocationStep"),
            workflow_step_id: this.string(null).nullable(),

            //ORM only
            //workflow_step: this.belongsTo(StoredWorkflowStep, 'workflow_step_id') //TODO create model
            //invocation: this.belongsTo(WorkflowInvocation, '') //TODO no backreference
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: ':url'
        },
        methods: {
            //TODO Galaxy api currently doesn't support steps fetch
            //$fetch: {
            //    name: 'fetch',
            //    http: {
            //        url: '',
            //        method: 'get',
            //    },
            //},
            $get: {
                name: 'get',
                http: {
                    url: '/:id',
                    method: 'get',
                },
            },
            //$create: {
            //    name: 'create',
            //    http: {
            //        url: '',
            //        method: 'post',
            //    },
            //},
            $update: {
                name: 'update',
                http: {
                    url: '/:id',
                    method: 'put',
                },
            },
            //$delete: {
            //    name: 'delete',
            //    http: {
            //        url: '/:id',
            //        method: 'delete',
            //    },
            //},
        }
    }
}

class WorkflowInvocation extends Common.Model {
    static entity = 'WorkflowInvocation';
    static primaryKey = 'id';
    static end_states = ['scheduled','error'];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            update_time: this.string(null).nullable(),
            uuid: this.string(null).nullable(),
            outputs: this.attr({}), //this.hasMany(HistoryDatasetAssociation, 'id'), //TODO
            output_collections: this.attr({}), //this.hasMany(HistoryDatasetCollectionAssociation, 'id'), //TODO
            history_id: this.string(null).nullable(),
            workflow_id: this.string(null).nullable(),
            state: this.string(null).nullable(),
            model_class: this.string("WorkflowInvocation"),
            inputs: this.attr({}),
            steps: this.hasMany(WorkflowInvocationStep, 'id'),

            //ORM only
            workflow: this.belongsTo(StoredWorkflow, 'workflow_id'),
            history: this.belongsTo(History, 'history_id'),
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: ':url/invocations',
            onResponse(response) {
                //TODO Bandaid to fix incorrect workflow_id
                let id = response.config.url.match(/\/api\/workflows\/([^/]+)\/invocations/);
                if (id) {
                    let data = response.data;
                    if (!(data instanceof Array)) data = [data];
                    for (let datum of data) {
                        datum.workflow_id = id[1];
                    }
                }
                return response.data;
            }
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
                    url: '',
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
    };
}

class StoredWorkflow extends Common.Model {
    static entity = 'StoredWorkflow';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            name: this.string("Unnamed workflow"),
            tags: this.attr([]),
            deleted: this.boolean(false),
            latest_workflow_uuid: this.string(null).nullable(),
            show_in_tool_panel: this.boolean(false),
            url: this.string(null).nullable(),
            number_of_steps: this.number(0),
            published: this.boolean(false),
            owner: this.string(null).nullable(),
            model_class: this.string("StoredWorkflow"),
            inputs: this.attr({}),
            annotation: this.string(''),
            version: this.number(0),
            steps: this.attr({}),

            //ORM only
            invocations: this.hasMany(WorkflowInvocation, 'workflow_id'),
        }
    }

    //TODO workflow menu logic
    //TODO GET /api/workflows/{encoded_workflow_id}/versions
    //TODO GET /api/workflows/{encoded_workflow_id}/download
    //TODO POST /api/workflows/import

    //TODO POST /api/workflows/{encoded_workflow_id}/invocations

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/workflows',
            //onResponse(response) {
            //    //TODO Bandaid to fix fetching shared workflows
            //    let found = response.config.url.match(/\/api\/workflows\/menu/);
            //    if (found) {
            //        return response.data.workflows;
            //    }
            //    return response.data;
            //}
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
                    url: '',
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
    database.register(WorkflowInvocationStep, Module);
    database.register(WorkflowInvocation, Module);
    database.register(StoredWorkflow, Module);
}

export {
    Module,
    StoredWorkflow,
    WorkflowInvocation,
    WorkflowInvocationStep,
    register,
};