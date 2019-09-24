import * as Common from "./_common";
import axios from "axios";


class Users extends Common.Model {
    static entity = 'users';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            username: this.string(null).nullable(),
            quota_percent: this.number(null).nullable(),
            preferences: this.attr({}),
            total_disk_usage: this.number(0),
            deleted: this.boolean(false),
            purged: this.boolean(false),
            nice_total_disk_usage: this.string("None"),
            quota: this.string("unlimited"),
            email: this.string(null).nullable(),
            is_admin: this.boolean(false),
            tags_used: this.attr([]),
        }
    }

    static async getCurrent() {
        let response = await this.$get({
            params: {
                id: 'current',
            }
        });
        return this.find(response.id);
    }

    static async registerUser(username, password, email) {
        // This is a total hack until user registration is enabled via the API
        const htmlconf = {...this.methodConf.http, responseType: 'text'};
        let response = await axios.get('/root/login', htmlconf);
        let csrf = response.data.match(/"session_csrf_token": "(\w+)"/);
        if (!csrf || !csrf[1]) throw Error('CSRF token could not be found');
        csrf = csrf[1];

        response = await axios.post('/user/create', {
            "disableCreate": true,
            "email": email,
            "password": password,
            "username": username,
            "confirm": password,
            "subscribe": null,
            "messageText": "",
            "messageVariant": "danger",
            "session_csrf_token": csrf,
            "isAdmin": false
        }, this.methodConf.http);
    }

    //TODO add methods to operate on rest of users api
    //GET /api/users/deleted Displays a collection (list) of deleted users.
    //GET /api/users/deleted/{encoded_id}
    //GET /api/users/{id}/information/inputs
    //PUT /api/users/{id}/information/inputs
    //GET /api/users/{id}/custom_builds
    //PUT /api/users/{id}/custom_builds/{key}
    //DELETE /api/users/{id}/custom_builds/{key}

    //GET /api/authenticate/baseauth
    static async getAPIKey(username, password, method = 'baseauth') {
        let response;
        switch (method) {
            case 'baseauth':
                response = await axios.get('/api/authenticate/baseauth', {
                    ...this.methodConf.http,
                    auth: {username, password}
                });
                return response.data.api_key;
            default:
                throw Error(method + ' authentication method not implemented');
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/users'
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
    database.register(Users, Module);
}

export {
    Users,
    Module,
    register,
};