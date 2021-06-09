/*
Configuration options specific to the app
 */

import {publicPath} from "../vue.config";

export const workflow_tag = "islandcompare";
export const workflow_owner = "brinkmanlab";
export const application_tag = "IslandCompare";
export const base_path = publicPath;
export let galaxy_path = `${window.location.protocol}//galaxy.islandcompare.ca`; // Default backend

if (process.env.NODE_ENV === 'production') {
    // Production only
    const host_parts = window.location.hostname.split('.');
    if (host_parts.length === 2) host_parts.unshift('galaxy'); // Handle root domain
    else host_parts[0] = 'galaxy'; // Replace subdomain otherwise
    galaxy_path = `${window.location.protocol}//${host_parts.join('.')}`;
} else {
    // Dev only
    //galaxy_path = `${window.location.protocol}//galaxy.dev.islandcompare.ca`;
    //galaxy_path = `${window.location.protocol}//localhost:8000`;
    //galaxy_path = `${window.location.protocol}//galaxy.stage.islandcompare.ca`;
    galaxy_path = `${window.location.protocol}//galaxy.aws.islandcompare.ca`;
}
