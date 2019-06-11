/*

 */
let uuidPromise = null;
let fetchedHistories = null;
let fetchedWorkflows = null;
let fetchedInvocations = null;

import { galaxy_load } from "@/store";
import { workflow_name } from "@/app.config";
import { getOrCreateUUID } from "@/auth";

export async function buildORM() {
    let galaxy = await galaxy_load;
    if (uuidPromise === null) uuidPromise = getOrCreateUUID();
    await uuidPromise;
    //Only fetch once across entire application for lifetime of window
    if (fetchedHistories === null)
        fetchedHistories = galaxy.histories.History.$fetch();
    if (fetchedWorkflows === null)
        fetchedWorkflows = galaxy.workflows.StoredWorkflow.$fetch({query: {show_published: true}});
    return galaxy;
}

export async function getConfiguredWorkflow() {
    // Load the workflow and all its components
    let galaxy = await galaxy_load;
    await buildORM();
    await fetchedWorkflows;
    let workflow = galaxy.workflows.StoredWorkflow.query().where('name', workflow_name).first();
    if (!workflow) {
        throw "IslandCompare workflow could not be found";
    }

    // Fetch each invocation individually by history id
    await fetchedHistories;
    if (fetchedInvocations === null) fetchedInvocations = Promise.all(galaxy.histories.History.all().map(h=>{
        if (h.tags.includes(workflow.id)) {
            return galaxy.workflows.WorkflowInvocation.$fetch({
                params: {url: workflow.url},
                query: {view: "element", step_details: true, history_id: h.id}
            });
        } else {
            return Promise.resolve();
        }
    }));

    await fetchedInvocations;

    let result = galaxy.workflows.StoredWorkflow.query().with('invocations', invocations => { //TODO break up this query across relevant components
        invocations.whereHas('history', history => {
            history.where('deleted', false).where('tags', tags => tags.includes(workflow.id));
        }).with('workflow').with('history');
    }).find(workflow.id);
    return result || workflow; //TODO clean up logic above to allow for empty workflows
}

export async function getUploadHistory() {
    // Load the user_data history and all its datasets
    let galaxy = await galaxy_load;
    await buildORM();
    await fetchedHistories;
    let history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).first();
    if (!history) {
        let response = await galaxy.histories.History.$create({
            data: {
                name: "Uploaded data", //TODO replace with app name
            }
        });
        history = galaxy.histories.History.find(response.id);
        history.tags.push('user_data');
        history.upload();
    } else {
        await galaxy.history_contents.HistoryDatasetAssociation.$fetch({
            params: {
                url: history.contents_url,
            }
        });
        await galaxy.history_contents.HistoryDatasetCollectionAssociation.$fetch({
            params: {
                url: history.contents_url,
            }
        });
    }
    return history;
}

export async function invokeConfiguredWorkflow(datasets, label, params) {
    // Invoke workflow, creating a new history and adding selected datasets to a collection first
    // All input is expected to be valid
    let galaxy = await galaxy_load;
    let workflow = await getConfiguredWorkflow();

    let response = null;
    //Create history to store run
    try {
        response = await galaxy.histories.History.$create({
            data: {
                name: label,
            }
        });
    } catch (e) {
        throw "Failed to create job history.";
    }

    let run_history = galaxy.histories.History.find(response.id);
    if (!run_history) {
        throw "Failed to create a job history.";
    }
    run_history.tags.push(workflow.id);
    run_history.upload();

    //Create collection of inputs in new history
    try {
        response = await galaxy.history_contents.HistoryDatasetCollectionAssociation.$create({
            params: {
                url: run_history.contents_url,
            },
            data: {
                name: "Selected datasets",
                type: 'dataset_collection',
                collection_type: 'list',
                //copy_elements: true, //TODO uncomment and test this, users can delete datasets during job run
                element_identifiers: datasets.map(model => ({
                    src: (model instanceof galaxy.history_contents.HistoryDatasetAssociation) ? model.hda_ldda : 'hdca', //TODO else 'hdca' is fragile
                    name: model.name + '_' + model.hid,
                    id: model.id,
                })),
            }
        });
    } catch (e) {
        run_history.delete();
        throw "Failed to create job dataset collection.";
    }
    //let params = Object.entries(this.params).reduce((a,[k,v])=>{a[k]=(v.toString ? v.toString() : v); return a}, {}); //https://github.com/galaxyproject/galaxy/issues/7654
    //Invoke workflow //TODO refactor this into the StoredWorkflow model
    try {
        response = await galaxy.workflows.WorkflowInvocation.$create({
            params: {
                url: workflow.url,
            },
            data: {
                //workflow_id: this.workflow.id,
                history_id: run_history.id,
                parameters: {
                    //TODO dynamically generate simple inputs, currently hardcoded to islandcompare
                    0:{"input":{"values":[{"src":"hdca","tags":[],"hid":8,"id":response.id}],"batch":false}},
                    13: {
                        cond: `c5-c4>${params.minimum_island_size}`,
                        "header_lines": "0"
                    },
                    27: {
                        "envs_0|name": "minimum_homologous_region",
                        "envs_0|val": params.minimum_homologous_region.toString(),
                        "envs_1|name": "min_cluster_size",
                        "envs_1|val": params.minimum_cluster_size.toString(),
                    }
                },
                parameters_normalized: true,
                batch: true,
                no_add_to_history: true,
            }
        });
        return galaxy.workflows.WorkflowInvocation.find(response.id);
    } catch (e) {
        run_history.delete();
        throw "Failed to create job.";
    }
}