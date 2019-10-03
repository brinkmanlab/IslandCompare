/*
User identity resolution code
 */
import Vue from 'vue'
import uuidgen from 'uuid/v1';

import { galaxy_load } from "@/store";
import { User } from "@/galaxy/src/api/users";

export let gidPromise = null;

export function getUUID() {
    let id = (new URLSearchParams(location.search)).get('uuid');
    if (!id) {
        let m = document.cookie.match(/galaxysession_user_uuid=([^;]+)/);
        if (m) {
            //If uuid cookie set, use that.
            id = m[1];
        }
    }
    if (id && gidPromise === null) gidPromise = setGlobalID(id);
    return id;
}

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
    document.cookie = `galaxysession_user_uuid=${id};path=/;max-age=31536000`;

    return id;
}

export async function setGlobalID(id) {
    // register filter to append ?uuid= to urls
    Vue.filter('auth', value=>value + (value.includes('?') ? '&' : '?') + 'key=' + id);

    // Set ?key= for all api requests
    const galaxy = await galaxy_load;
    galaxy.applyMethodConf(conf=>{
        if (conf.http.params) conf.http.params.key = id;
        else conf.http.params = {key: id};
    });

    return id;
}