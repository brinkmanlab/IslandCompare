import * as Common from "./_common";


class Genome extends Common.Model {
    static entity = 'Genome';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null),
            name: this.string(null).nullable(),
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/genomes',
            onResponse(response) {
                //TODO Bandaid to deal with lazy api output
                if (response.config.url.match('/api/genomes$')) {
                    let data = [];
                    for (let i = 0; i < response.data.length; ++i) {
                        if (response.data[i][0] === '----- Additional Species Are Below -----') continue; // Handle https://github.com/galaxyproject/galaxy/issues/8612
                        data[i] = {name: response.data[i][0], id: response.data[i][1]};
                    }
                    return data;
                }
                return response.data;
            },
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
            /*$create: {
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
            },*/
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
    database.register(Genome, Module);
}

export {
    Genome,
    Module,
    register,
};