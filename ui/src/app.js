/*
 * Code specific to IslandCompare app
 * Contains functions for initialising state and populating ORM
 */
let uuidPromise = null;
let fetchedHistories = null;
let fetchedWorkflows = null;
let fetchedInvocations = null;

import { api as galaxy } from 'galaxy-client'
import { workflow_name } from "@/app.config";
import { getOrCreateUUID } from "@/auth";

export async function buildORM() {
    if (uuidPromise === null) uuidPromise = getOrCreateUUID();
    await uuidPromise;
    //Only fetch once across entire application for lifetime of window
    if (fetchedHistories === null)
        fetchedHistories = galaxy.histories.History.fetch();
    if (fetchedWorkflows === null)
        fetchedWorkflows = galaxy.workflows.StoredWorkflow.fetch({params: {show_published: true}});
    return galaxy;
}

export async function getConfiguredWorkflow() {
    // Load the workflow and all its components
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
        await buildORM();
        await fetchedHistories;
        const workflow = await workflowPromise;
        const histories = galaxy.histories.History.query().where('deleted', false).where('tags',  tags=>tags.includes(workflow.id)).get();
        fetchedInvocations = workflow.fetch_invocations(histories);
    }
    return await fetchedInvocations;
}

export async function getUploadHistory() {
    // Load the user_data history and all its datasets
    await buildORM();
    await fetchedHistories;
    let history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).with('datasets.history').first();
    if (!history) {
        history = await galaxy.histories.History.post({
            name: "Uploaded data", //TODO replace with app name
        });
        history.tags.push('user_data');
        history.put(['tags']);
    } else {
        history.loadContents();
        //history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).first(); //TODO .with('datasets.history').with('collections.history')
    }
    return history;
}
