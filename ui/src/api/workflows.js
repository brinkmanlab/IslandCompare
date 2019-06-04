import * as Common from "./_common";
//import axios from "axios";
import {HistoryDatasetAssociation, HistoryDatasetCollectionAssociation } from "./history_contents";
import { History } from "./histories";
import { Job } from "./jobs";


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
            jobs: this.hasMany(Job, 'workflow_invocation_step_id'),
            job_id: this.string(null).nullable(),
            outputs: this.attr(), //this.hasMany(HistoryDatasetAssociation, 'id'), //TODO
            order_index: this.number(0),
            output_collections: this.attr(), //this.hasMany(HistoryDatasetCollectionAssociation, 'id'), //TODO
            workflow_step_label: this.string(''),
            state: this.string(null).nullable(),
            action: this.string(null).nullable(),
            model_class: this.string("WorkflowInvocationStep"),
            workflow_step_id: this.string(null).nullable(),
            workflow_invocation_id: this.string(null),

            //ORM only
            //workflow_step: this.belongsTo(StoredWorkflowStep, 'workflow_step_id') //TODO create model
            invocation: this.belongsTo(WorkflowInvocation, 'workflow_invocation_id') //TODO no backreference
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

/*class WorkflowInvocationOutput extends Common.Model {
    static entity = "WorkflowInvocationOutput";
    static primaryKey = 'label';

    static fields() {
        return {
            label: this.string(null),
            id: this.string(null),
            src: this.string(null),
        }
    }
}*/

class WorkflowInvocation extends Common.Model {
    static entity = 'WorkflowInvocation';
    static primaryKey = 'id';
    static end_states = ["cancelled", "error", "done"];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            update_time: this.string(null).nullable(),
            uuid: this.string(null).nullable(),
            outputs: this.attr({}),
            output_collections: this.attr({}),
            history_id: this.string(null).nullable(),
            workflow_id: this.string(null).nullable(),
            state: this.string(null).nullable(),
            model_class: this.string("WorkflowInvocation"),
            inputs: this.attr({}),
            steps: this.hasMany(WorkflowInvocationStep, 'workflow_invocation_id'),

            //ORM only
            workflow: this.belongsTo(StoredWorkflow, 'workflow_id'),
            history: this.belongsTo(History, 'history_id'),
            output_models: this.hasManyBy(HistoryDatasetAssociation, 'outputs'),
            output_collection_models: this.hasManyBy(HistoryDatasetCollectionAssociation, 'output_collections'),
        }
    }

    get_base_url() {
        let workflow = this.workflow;
        if (!workflow) workflow = StoredWorkflow.find(this.workflow_id);
        return workflow.url;
    }

    states() {
        return this.steps.reduce((acc, cur) => {
            acc[cur.state] = (acc[cur.state] || 0) + 1;
            return acc
        }, {});
    }

    aggregate_state() {
        if (this.state === "cancelled") return this.state;
        let states = this.states();
        if (Object.entries(states).length === 0) return "new";
        if (states.error) return "error";
        if (states.new) return "running";
        return "done";
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

                //TODO Bandaid to deal with invocation not storing ids of steps
                if (Array.isArray(response.data)) {
                    response.data.forEach(invocation =>{
                        if (invocation.hasOwnProperty('steps')) invocation.steps.forEach(step => {
                        step.workflow_invocation_id = response.data.id
                    })});
                } else if (response.data.hasOwnProperty('steps')) {
                    response.data.steps.forEach(step => {
                        step.workflow_invocation_id = response.data.id;
                        if (step.jobs) {
                            step.jobs.forEach(job => {
                                job.workflow_invocation_step_id = step.id;
                            });
                        }
                    });
                }

                //TODO Bandaid to deal with receiving dict for outputs and output_collections
                //response.data.outputs = Object.keys(response.data.outputs).reduce((acc, cur)=>{acc.append({label: cur, ...response.data.outputs[cur]})}, []);
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