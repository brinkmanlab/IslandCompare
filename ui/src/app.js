/*
 * Code specific to IslandCompare app
 * Contains functions for initialising state and populating ORM
 */
let stateFetched = false; // Prevent fetching more than once
let historyFetched = false; // Prevent fetching more than once
export let invocationsFetched = false;

import { api as galaxy } from 'galaxy-client'
import {workflow_tag, application_tag, workflow_owner} from "./app.config";
import {getOrCreateUUID, getAPIKey, setGlobalKey} from "./auth";

// Fetch all required state from api
export async function fetchState(createUUID = false) {
    if (stateFetched) return;  //Only fetch once across entire application for lifetime of window
    stateFetched = true;
    let key = '';
    if (createUUID) key= await getOrCreateUUID();
    else key = await getAPIKey();
    setGlobalKey(key);
    if (!key) {
        stateFetched = false;
        throw "No api key";
    }

    await galaxy.histories.History.fetch();
    await galaxy.workflows.StoredWorkflow.fetch({params: {show_published: true}});

    // Workflow inputs are not loaded in fetch, load them.
    const workflow = getConfiguredWorkflow();
    if (!workflow) throw("Failed to load workflow");
    if (Object.keys(workflow.inputs).length === 0) {
        workflow.reload(); //Get input details
    }

    // Fetch workflow invocations for histories that are not deleted
    const histories = galaxy.histories.History.query().where('deleted', false).where('tags', tags => tags.includes(workflow.id)).get();
    workflow.fetch_invocations(histories).then(() => invocationsFetched = true);
}

export async function fetchStateAndUploadHistory(createUUID = false) {
    await fetchState(createUUID);

    // Load or create the user_data history and all its datasets
    let history = getUploadHistory();
    if (!history) {
        historyFetched = true;
        history = await galaxy.histories.History.post({
            name: "Uploaded data", //TODO replace with app name
        });
        history.tags.push('user_data');
        history.put(['tags']);
    } else if (history && !historyFetched) {
        historyFetched = true;
        history.loadContents();
    }
    return history;
}

export function getConfiguredWorkflow() {
    // Look for workflow with owner but fail back to any owner
    return galaxy.workflows.StoredWorkflow.query().where('owner', workflow_owner).where('tags', tags=>tags.includes(workflow_tag)).with('inputs|steps').first() ||
        galaxy.workflows.StoredWorkflow.query().where('tags', tags=>tags.includes(workflow_tag)).with('inputs|steps').first();
}

export function getInvocations(workflow) {
    return galaxy.workflows.WorkflowInvocation.query().whereHas('history', q => {
        q.where('deleted', false)
            .where('tags', tags=>tags.includes(application_tag) || tags.includes(workflow.id))
    }).with('history|workflow|steps.jobs').get();
}

export function getUploadHistory() {
    return galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).with('datasets.history').first();
}

export async function onInvocation(invocation) {
    invocation = await invocation;
    const history = galaxy.histories.History.find(invocation.history_id);
    history.tags.push(application_tag);
    history.put(['tags']);
}
