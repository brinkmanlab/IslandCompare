/*
User identity resolution code
 */
import Vue from 'vue'
import uuidgen from 'uuid/v1';
import axios from 'axios';
import { galaxy_path } from "@/app.config";

export function getUUID() {
    let id = (new URLSearchParams(location.search)).get('uuid');
    if (!id) {
        let m = document.cookie.match(/galaxysession_user_uuid=([^;]+)/);
        if (m) {
            //If uuid cookie set, use that.
            id = m[1];
        }
    }
    if (id) setGlobalID(id);
    return id;
}

export async function getOrCreateUUID() {
    let id = getUUID();
    if (!id) {
        //generate uuid
        // TODO move this logic to backend
        id = uuidgen();
        let tag = location.search.lastIndexOf('#');
        let url = "";
        if (tag >= 0) {
            url = location.search.slice(0, tag) + (location.search.includes('?') ? '&' : '?') + id + location.search.slice(tag);
        } else {
            url = location.search + (location.search.includes('?') ? '&' : '?') + "uuid=" + id;
        }
        history.replaceState(history.state, "Analysis", url);
        await setGlobalID(id);
        alert("Be sure to bookmark this page to return to your work. The URL is unique to you."); //TODO replace with a html popup
    }
    document.cookie = `galaxysession_user_uuid=${id};path=/;max-age=31536000`;

    /*/ fetch api key
    let response = await axios.get('/api/users', { //eslint-disable-line
        params: {
            uuid: id,
        }
    });

    response = await axios.get(`/api/users/${response.data[0].id}/api_key/inputs`);
    */



    return id;
}

export async function setGlobalID(id) {
    // register filter to append ?uuid= to urls
    Vue.filter('auth', value=>value + (value.includes('?') ? '&' : '?') + 'uuid=' + id);
    //this is a bandaid to get a session key from the galaxy frontend rather than an api key, api keys are not available to remote auth users
    let response = await axios.get('/user', { //eslint-disable-line
        baseURL: galaxy_path,
        params: {
            uuid: id,
        }
    });
}