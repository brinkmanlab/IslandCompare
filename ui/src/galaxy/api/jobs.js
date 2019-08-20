import * as Common from "./_common";
import {History} from "./histories";
import {HistoryDatasetAssociation} from "./history_contents";


class Job extends Common.Model {
    static entity = 'Jobs';
    static primaryKey = 'id';
    static end_states = ['ok', 'error', 'deleted', 'deleted_new'];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null),
            tool_id: this.string().nullable(),
            update_time: this.string().nullable(),
            history_id: this.string().nullable(),
            exit_code: this.number().nullable(),
            state:	this.string().nullable(), // ‘new’, ‘upload’, ‘waiting’, ‘queued’, ‘running’, ‘ok’, ‘error’, ‘paused’, ‘deleted’, ‘deleted_new’
            create_time: this.string().nullable(),
            model_class: this.string("Job"),
            inputs: this.attr(),
            outputs: this.attr(),
            params: this.attr(),

            //ORM Only
            history: this.belongsTo(History, 'history_id'),
            workflow_invocation_step_id: this.string().nullable(),
        }
    }

    async get_error_log(label) {
        let self = this;
        if (self.state !== 'error') return '';
        if (!self.inputs || !self.outputs) {
            await self.reload();
            self = self.constructor.find(self.id);
        }
        if (!self.history) {
            self.history = await History.findOrLoad(self.history_id);
        }
        let log = '';

        // Resolve input identifiers
        let input_identifier = Object.keys(self.inputs)
                .filter(key => self.params.hasOwnProperty(`${key}|__identifier__`))
                .map(key => self.params[`${key}|__identifier__`]);
        if (input_identifier.length === 1) input_identifier = input_identifier[0];
        else if (input_identifier.length > 1) input_identifier = `[${input_identifier.join(', ')}]`;
        else input_identifier = '';

        for (const [key, val] of Object.entries(self.outputs)) {
            switch (val.src) {
                case 'hda': {
                    const hda = await HistoryDatasetAssociation.findOrLoad(val.id, self.history.contents_url);
                    if (hda.state === 'error') log += `${label || self.tool_id} on ${input_identifier} - ${key}: ${hda.misc_info}\n`;
                    break;
                }
                case 'hdca': {
                    // TODO
                    break;
                }
                default:
                    break;
            }
        }
        return log;
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/jobs'
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
            //$create: {
            //    name: 'create',
            //    http: {
            //        url: '', //TODO
            //        method: 'post',
            //    },
            //},
            //$update: {
            //    name: 'update',
            //    http: {
            //        url: '/:id', //TODO
            //        method: 'put',
            //    },
            //},
            //$delete: {
            //    name: 'delete',
            //    http: {
            //        url: '/:id', //TODO
            //        method: 'delete',
            //    },
            //},
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
    database.register(Job, Module);
}

export {
    Job,
    Module,
    register,
};