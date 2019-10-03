/*
 * Code specific to IslandCompare app
 * Contains functions for initialising state and populating ORM
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
    let workflow = galaxy.workflows.StoredWorkflow.query().where('name', workflow_name).with('inputs|steps').first();
    if (!workflow) {
        throw Error(workflow_name + " workflow could not be found");
    }

    if (Object.keys(workflow.inputs).length === 0) {
        await workflow.reload(); //Get input details
        workflow = galaxy.workflows.StoredWorkflow.query().where('name', workflow_name).with('inputs|steps').first();
    }
    return workflow;
}

export async function getInvocations(workflowPromise) {
    if (fetchedInvocations === null) {
        let galaxy = await galaxy_load;
        await fetchedHistories;
        const workflow = await workflowPromise;
        const histories = galaxy.histories.History.all().filter(h => !h.deleted && h.tags.includes(workflow.id));
        fetchedInvocations = workflow.fetch_invocations(histories);
    }
    return await fetchedInvocations;
}

export async function getUploadHistory() {
    // Load the user_data history and all its datasets
    let galaxy = await galaxy_load;
    await buildORM();
    await fetchedHistories;
    let history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).with('datasets.history').first();
    if (!history) {
        let response = await galaxy.histories.History.$create({
            data: {
                name: "Uploaded data", //TODO replace with app name
            }
        });
        history = galaxy.histories.History.find(response.id);
        history.tags.push('user_data');
        history.post();
    } else {
        // TODO should this be moved to the model?
        if (history.datasets.length === 0) {
            await galaxy.history_contents.HistoryDatasetAssociation.$fetch({
                params: {
                    url: history.contents_url,
                }
            });
        }
        /*if (history.collections.length === 0) { TODO
            await galaxy.history_contents.HistoryDatasetCollectionAssociation.$fetch({
                params: {
                    url: history.contents_url,
                }
            });
        }*/
        //history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).first(); //TODO .with('datasets.history').with('collections.history')
    }
    return history;
}