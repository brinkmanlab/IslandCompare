import VuexORM, { Database } from '@vuex-orm/core'
import VuexORMAxios from "@vuex-orm/plugin-axios/src"; // TODO remove /src

const database = new Database();

//import * as annotations from './api/annotations';
//annotations.register(database);

//import * as cloud from './api/cloud';
//cloud.register(database);

//import * as configuration from './api/configuration';
//configuration.register(database);

//import * as dataset_collections from './api/dataset_collections';
//dataset_collections.register(database);

//import * as datasets from './api/datasets';
//datasets.register(database);

//import * as datatypes from './api/datatypes';
//datatypes.register(database);

//import * as extended_metadata from './api/extended_metadata';
//extended_metadata.register(database);

//import * as folder_contents from './api/folder_contents';
//folder_contents.register(database);

//import * as folders from './api/folders';
//folders.register(database);

//import * as forms from './api/forms';
//forms.register(database);

import * as genomes from './api/genomes';
genomes.register(database);

//import * as group_roles from './api/group_roles';
//group_roles.register(database);

//import * as group_users from './api/group_users';
//group_users.register(database);

//import * as groups from './api/groups';
//groups.register(database);

import * as histories from './api/histories';
histories.register(database);

import * as history_contents from './api/history_contents';
history_contents.register(database);

//import * as item_tags from './api/item_tags';
//item_tags.register(database);

//import * as job_files from './api/job_files';
//job_files.register(database);

import * as jobs from './api/jobs';
jobs.register(database);

//import * as libraries from './api/libraries';
//libraries.register(database);

//import * as library_contents from './api/library_contents';
//library_contents.register(database);

//import * as library_datasets from './api/library_datasets';
//library_datasets.register(database);

//import * as metrics from './api/metrics';
//metrics.register(database);

//import * as page_revisions from './api/page_revisions';
//page_revisions.register(database);

//import * as pages from './api/pages';
//pages.register(database);

//import * as plugins from './api/plugins';
//plugins.register(database);

//import * as provenance from './api/provenance';
//provenance.register(database);

//import * as quotas from './api/quotas';
//quotas.register(database);

//import * as remote_files from './api/remote_files';
//remote_files.register(database);

//import * as roles from './api/roles';
//roles.register(database);

//import * as search from './api/search';
//search.register(database);

//import * as tool_data from './api/tool_data';
//tool_data.register(database);

//import * as tool_dependencies from './api/tool_dependencies';
//tool_dependencies.register(database);

//import * as tool_shed_repositories from './api/tool_shed_repositories';
//tool_shed_repositories.register(database);

//import * as tools from './api/tools';
//tools.register(database);

//import * as toolshed from './api/toolshed';
//toolshed.register(database);

//import * as tours from './api/tours';
//tours.register(database);

//import * as uploads from './api/uploads';
//uploads.register(database);

import * as users from './api/users';
users.register(database);

//import * as visualizations from './api/visualizations';
//visualizations.register(database);

//import * as webhooks from './api/webhooks';
//webhooks.register(database);

import * as workflows from './api/workflows';
workflows.register(database);

import * as common from './api/_common';

// Call this function after async import() of this file. This file is expensive to load.
//TODO VuexORM.Database.register is syncrounous and can be called after load to register models on the fly.
//      database.start() needs to be called after each batch of models are registered.
function register(store, http={}) {
    VuexORM.use(VuexORMAxios, {
        database,
        http: {
            baseURL: '/',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json',
            },
            ...http,
        },
    });
    VuexORM.install(database, {namespace: 'galaxy'})(store)
}

/**
 * Apply a function to all model methodConf properties
 * Use this to update global axios settings, like the api key
 * @param func function that takes a single argument, the methodConf instance of the model
 */
function applyMethodConf(func) {
    for (const submodule of Object.values(require('.'))) {
        for (const model of Object.values(submodule)) {
            if (model.hasOwnProperty('methodConf')) {
                func(model.methodConf);
            }
        }
    }
}

export {
    database,
    register,
    applyMethodConf,
    common,
//    annotations,
//    authenticate,
//    cloud,
//    configuration,
//    dataset_collections,
//    datasets,
//    datatypes,
//    extended_metadata,
//    folder_contents,
//    folders,
//    forms,
    genomes,
//    group_roles,
//    group_users,
//    groups,
    histories,
    history_contents,
//    item_tags,
//    job_files,
    jobs,
//    libraries,
//    library_contents,
//    library_datasets,
//    metrics,
//    page_revisions,
//    pages,
//    plugins,
//    provenance,
//    quotas,
//    remote_files,
//    roles,
//    search,
//    tool_data,
//    tool_dependencies,
//    tool_shed_repositories,
//    tools,
//    toolshed,
//    tours,
//    uploads,
    users,
//    visualizations,
//    webhooks,
    workflows,
};