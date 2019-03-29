import * as Common from "./_common";

class ToolSection extends Common.Model {
    static entity = 'ToolSection';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null),
            model_class: this.string("ToolSection"),
            version: this.string(null).nullable(),
            name: this.string("Unnamed"),
            elems: this.hasMany(Tool, 'id'),
        }
    }

    //TODO API doesn't provide section endpoints
    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/tools'
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '',
                    method: 'get',
                },
            },
            /*$get: {
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
            },*/
        }
    }
}

class Tool extends Common.Model {
    static entity = 'Tool';
    static primaryKey = 'id';
    static parent = 'section';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null),
            panel_section_name: this.string(null).nullable(),
            config_file: this.string(null).nullable(),
            description: this.string(null).nullable(),
            labels: this.attr([]),
            edam_operations: this.attr([]),
            form_style: this.string("regular"),
            edam_topics: this.attr([]),
            panel_section_id: this.string(null).nullable(),
            version: this.string(null).nullable(),
            link: this.string(null).nullable(),
            target: this.string(null).nullable(),
            min_width: this.number(-1),
            model_class: this.string("Tool"),
            name: this.string(null).nullable(),

            //ORM Only
            section: this.belongsTo(ToolSection, 'panel_section_id')
        }
    }

    //TODO GET /api/tools/{tool_id}/build
    //TODO GET /api/tools/{tool_id}/test_data_path?tool_version={tool_version}
    //TODO GET /api/tools/{tool_id}/test_data_download?tool_version={tool_version}&filename={filename}
    //TODO GET /api/tools/tests_summary
    //TODO GET /api/tools/{tool_id}/test_data?tool_version={tool_version}
    //TODO GET /api/tools/{tool_id}/reload
    //TODO GET /api/tools/all_requirements
    //TODO GET /api/tools/{tool_id}/requirements
    //TODO POST /api/tools/{tool_id}/dependencies
    //TODO DELETE /api/tools/{tool_id}/dependencies
    //TODO POST /api/tools/{tool_id}/build_dependency_cache
    //TODO GET /api/tools/{tool_id}/diagnostics
    //TODO GET /api/tools/error_stack

    //TODO api doesn't provide method to fetch tools independent of section
    static $fetch(options) {
        return ToolSection.$fetch(options);
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/tools'
        },
        methods: {
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
            $create: {
                name: 'create',
                http: {
                    url: '',
                    method: 'post',
                },
            },
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
    database.register(Tool, Module);
    database.register(ToolSection, Module);
}

export {
    Tool,
    ToolSection,
    Module,
    register,
};