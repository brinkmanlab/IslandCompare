/*
User identity resolution code
 */
import Vue from 'vue'
import uuidgen from 'uuid/v1';

import {api, config} from 'galaxy-client'
const User = api.users.User;

export let gidPromise = null; // Only set uuid globally once per page load

/**
 * Attempt to recover uuid from url or browser cache
 * @returns {string} UUID of user
 */
export async function getUUID() {
    let id = (new URLSearchParams(location.search)).get('uuid') || localStorage.getItem('galaxysession_user_uuid');
    if (id && !await User.verifyAPIKey(id)) {
        // id did not validate against server
        id = '';
    }
    if (id) {
        if (gidPromise === null) gidPromise = setGlobalID(id);
        await gidPromise;
    }
    return id;
}

/**
 * Helper to insert the uuid into the currentl url query
 * @param router Router instance of the Vue instance
 * @param route Route instance of the Vue instance
 * @returns {Promise<void>}
 */
export async function updateRoute(router, route) {
    // Force uuid into url when navigating to this page
    const uuid = await getUUID();
    if (!('uuid' in route.query) && uuid) {
        router.replace({query: {uuid: uuid, ...route.query}}).catch(err=>err); // Triggers a duplicate navigation exception, likely doesn't check the query
    }
}

/**
 * Generate a UUID for the user if it does not already exist
 * @returns {Promise<string>} UUID of user
 */
export async function getOrCreateUUID() {
    let id = await getUUID();
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
        history.replaceState(history.state, "Analysis", `${location.origin}${location.pathname}?${search}${location.hash}`);
    }

    // Set ?key= for all api requests
    if (config.params) config.params.key = id;
    else config.params = {key: id};

    return id;
}