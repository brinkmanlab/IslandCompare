import VuexORM, { Database } from '@vuex-orm/core'
import VuexORMAxios from '@vuex-orm/plugin-axios'

const database = new Database();

import * as annotations from './api/annotations';
annotations.register(database);

import * as authenticate from './api/authenticate';
authenticate.register(database);

import * as cloud from './api/cloud';
cloud.register(database);

import * as configuration from './api/configuration';
configuration.register(database);

import * as dataset_collections from './api/dataset_collections';
dataset_collections.register(database);

import * as datasets from './api/datasets';
datasets.register(database);

import * as datatypes from './api/datatypes';
datatypes.register(database);

import * as extended_metadata from './api/extended_metadata';
extended_metadata.register(database);

import * as folder_contents from './api/folder_contents';
folder_contents.register(database);

import * as folders from './api/folders';
folders.register(database);

import * as forms from './api/forms';
forms.register(database);

import * as genomes from './api/genomes';
genomes.register(database);

import * as group_roles from './api/group_roles';
group_roles.register(database);

import * as group_users from './api/group_users';
group_users.register(database);

import * as groups from './api/groups';
groups.register(database);

import * as histories from './api/histories';
histories.register(database);

import * as history_contents from './api/history_contents';
history_contents.register(database);

import * as item_tags from './api/item_tags';
item_tags.register(database);

import * as job_files from './api/job_files';
job_files.register(database);

import * as jobs from './api/jobs';
jobs.register(database);

import * as libraries from './api/libraries';
libraries.register(database);

import * as library_contents from './api/library_contents';
library_contents.register(database);

import * as library_datasets from './api/library_datasets';
library_datasets.register(database);

import * as metrics from './api/metrics';
metrics.register(database);

import * as page_revisions from './api/page_revisions';
page_revisions.register(database);

import * as pages from './api/pages';
pages.register(database);

import * as plugins from './api/plugins';
plugins.register(database);

import * as provenance from './api/provenance';
provenance.register(database);

import * as quotas from './api/quotas';
quotas.register(database);

import * as remote_files from './api/remote_files';
remote_files.register(database);

import * as roles from './api/roles';
roles.register(database);

import * as search from './api/search';
search.register(database);

import * as tool_data from './api/tool_data';
tool_data.register(database);

import * as tool_dependencies from './api/tool_dependencies';
tool_dependencies.register(database);

import * as tool_shed_repositories from './api/tool_shed_repositories';
tool_shed_repositories.register(database);

import * as tools from './api/tools';
tools.register(database);

import * as toolshed from './api/toolshed';
toolshed.register(database);

import * as tours from './api/tours';
tours.register(database);

import * as uploads from './api/uploads';
uploads.register(database);

import * as users from './api/users';
users.register(database);

import * as visualizations from './api/visualizations';
visualizations.register(database);

import * as webhooks from './api/webhooks';
webhooks.register(database);

import * as workflows from './api/workflows';
workflows.register(database);

import uuid from 'uuid/v1';
//import axios from 'axios';

async function getDatabase() {
    let id = (new URLSearchParams(location.search)).get('uuid');
    if (!id) {
        id = uuid();
        document.cookie = `galaxysession_user_uuid=${id};path=/;max-age=31536000`;
        let tag = location.search.lastIndexOf('#');
        if (tag >= 0) {
            location.search = location.search.slice(0, tag) + (location.search.includes('?') ? '&' : '?') + id + location.search.slice(tag);
        } else {
            location.search += (location.search.includes('?') ? '&' : '?') + "uuid=" + id;
        }
    }
    /*/ fetch api key
    let response = await axios.get('/api/users', { //eslint-disable-line
        params: {
            uuid: id,
        }
    });

    response = await axios.get(`/api/users/${response.data[0].id}/api_key/inputs`);
    */
    VuexORM.use(VuexORMAxios, {
        database,
        http: {
            baseURL: '/',
            params: {
                uuid: id,
                //key: 'admin',
            },
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json',
            }
        },
    });

    return database;
}

import * as common from './api/_common';

export {
    getDatabase,

    common,
    annotations,
    authenticate,
    cloud,
    configuration,
    dataset_collections,
    datasets,
    datatypes,
    extended_metadata,
    folder_contents,
    folders,
    forms,
    genomes,
    group_roles,
    group_users,
    groups,
    histories,
    history_contents,
    item_tags,
    job_files,
    jobs,
    libraries,
    library_contents,
    library_datasets,
    metrics,
    page_revisions,
    pages,
    plugins,
    provenance,
    quotas,
    remote_files,
    roles,
    search,
    tool_data,
    tool_dependencies,
    tool_shed_repositories,
    tools,
    toolshed,
    tours,
    uploads,
    users,
    visualizations,
    webhooks,
    workflows,
}