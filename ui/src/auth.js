/*
User identity resolution code
 */
import Vue from 'vue'
import uuidgen from 'uuid/v1';

import {applyMethodConf} from "@/galaxy/src";
import { User } from "@/galaxy/src/api/users";

export let gidPromise = null; // Only set uuid globally once per page load

/**
 * Attempt to recover uuid from url or browser cache
 * @returns {string} UUID of user
 */
export function getUUID() {
    let id = (new URLSearchParams(location.search)).get('uuid') || localStorage.getItem('galaxysession_user_uuid');
    if (id && gidPromise === null) gidPromise = setGlobalID(id);
    return id;
}

/**
 * Generate a UUID for the user if it does not already exist
 * @returns {Promise<string>} UUID of user
 */
export async function getOrCreateUUID() {
    let id = getUUID();
    if (!id) {
        //generate uuid and get api key
        const uuid = uuidgen();
        await User.registerUser(uuid, uuid, uuid + '@external.ex');
        id = await User.getAPIKey(uuid + '@external.ex', uuid);
        if (gidPromise === null) gidPromise = setGlobalID(id);
        await gidPromise;
        alert("Be sure to bookmark this page to return to your work. The URL is unique to you."); //TODO replace with a html popup
    }
    localStorage.setItem('galaxysession_user_uuid', id);

    return id;
}

/**
 * Helper to register the UUID for the user throughout the app
 * @param id {string} User UUID
 * @returns {Promise<string>} Forwards the passed in UUID value
 */
export async function setGlobalID(id) {
    // register filter to append ?uuid= to urls
    Vue.filter('auth', value=>value + (value.includes('?') ? '&' : '?') + 'key=' + id);

    const search = new URLSearchParams(location.search);
    if (!search.get('uuid')) {
        // Rewrite url to contain uuid
        search.append('uuid', id);
        history.replaceState(history.state, "Analysis", `${location.origin}${location.pathname}${search}${location.hash}`);
    }

    // Set ?key= for all api requests
    applyMethodConf(conf=>{
        if (conf.http.params) conf.http.params.key = id;
        else conf.http.params = {key: id};
    });

    return id;
}