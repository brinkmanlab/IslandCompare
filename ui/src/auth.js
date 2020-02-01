/*
User identity resolution code
 */
import Vue from 'vue'
import uuidgen from 'uuid/v1';

import {api, config} from 'galaxy-client'
const User = api.users.User;

// If you change any of these values you should add a conditional to getAPIKey to ensure that returning users can recover their account
export const api_key_key = 'galaxy_api_key';
export const user_id_key = 'user_id';
export const user_id_param = 'id';

/**
 * Attempt to recover API key from url or browser cache
 * @returns {string} API key of user
 */
export async function getAPIKey() {
    let id = (new URLSearchParams(location.search)).get(user_id_param);
    const stored_id = localStorage.getItem(user_id_key);
    let old = false;
    let key = '';
    // Allow overriding stored key with id in url param
    if (!id || id === stored_id) key = localStorage.getItem(api_key_key);
    if (!key) {
        // Key not stored, get using id
        id = id || stored_id;
        if (id) {
            key = await User.getAPIKey(id + '@external.ex', id); // TODO This needs to be changed to get username/password when authentication is added
            if (key) localStorage.setItem(api_key_key, key);
        }
    }
    if (!key) {
        // Check for old way of storing API key
        // TODO remove this after all the testers are done with their data
        key = (new URLSearchParams(location.search)).get('uuid') || localStorage.getItem('galaxysession_user_uuid');
        if (key) {
            old = true;
        }
    }
    if (key && !await User.verifyAPIKey(key)) {
        // key did not validate against server
        key = '';
        // Purge stored values
        localStorage.removeItem(user_id_key);
        localStorage.removeItem(api_key_key);
        localStorage.removeItem('galaxysession_user_uuid');
    }
    if (key) {
        if (id) {
            // Ensure id is stored
            localStorage.setItem(user_id_key, id);
        }
        if (old) {
            // Change to new storage method
            const user = await User.getCurrent({params:{key: key}}); // Provide api key as it is not globally registered
            console.log("Recovering user id: " + user.username); //eslint-disable-line
            if (user.username) {
                localStorage.setItem(user_id_key, user.username);
                localStorage.removeItem('galaxysession_user_uuid');
            }
        }
    }
    return key;
}

/**
 * Helper to insert the uuid into the current url query
 * @param router Router instance of the Vue instance
 * @param route Route instance of the Vue instance
 * @returns {Promise<void>}
 */
export async function updateRoute(router, route) {
    // Force uuid into url when navigating to this page
    const id = localStorage.getItem(user_id_key);
    if (!id) console.log("Failed to update route, no id found."); //eslint-disable-line
    if (!(user_id_param in route.query) && id) {
        router.replace({query: {[user_id_param]: id, ...route.query}}).catch(err=>err); // Triggers a duplicate navigation exception, likely because it doesn't check the query for changes
    }
}

/**
 * Generate a UUID for the user if it does not already exist and get API key
 * @returns {Promise<string>} API key of user
 */
export async function getOrCreateUUID() {
    let key = await getAPIKey();
    if (!key) {
        //generate uuid and get api key
        const uuid = uuidgen();
        localStorage.setItem(user_id_key, uuid);
        await User.registerUser(uuid, uuid, uuid + '@external.ex');
        key = await User.getAPIKey(uuid + '@external.ex', uuid);
        if (key) {
            alert("Be sure to bookmark this page to return to your work. The URL is unique to you."); //TODO replace with a html popup
        }
    }

    return key;
}

/**
 * Helper to register the API key for the user throughout the app
 * @param key {string} API key
 * @returns {string} API key
 */
export function setGlobalKey(key) {
    if (!key) return;
    // register filter to append ?uuid= to urls
    Vue.filter('auth', value=>value + (value.includes('?') ? '&' : '?') + 'key=' + key);

    localStorage.setItem(api_key_key, key);

    // Set ?key= for all api requests
    if (config.params) config.params.key = key;
    else config.params = {key: key};

    return key;
}

/*
const search = new URLSearchParams(location.search);
const uuid = search.get('uuid');
if (!uuid || uuid !== key) {
    // Rewrite url to contain uuid
    search.set('uuid', key);
    history.replaceState(history.state, "Analysis", `${location.origin}${location.pathname}?${search}${location.hash}`);
}
 */
