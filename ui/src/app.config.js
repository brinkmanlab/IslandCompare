/*
Configuration options specific to the app
 */

import {publicPath} from "../vue.config";

export const workflow_tag = "islandcompare";
export const workflow_owner = "brinkmanlab";
export const application_tag = "IslandCompare";
export const base_path = publicPath;
export let galaxy_path = `${window.location.protocol}//galaxy.islandcompare.ca`; // Default backend

const subdomain = /^islandcompare/;
const fulldomain = /^islandcompare.*\..*\..*$/;

if (process.env.NODE_ENV === 'production') {
    // Production only
    if (fulldomain.test(window.location.hostname)) {
        galaxy_path = `${window.location.protocol}//${window.location.hostname.replace(subdomain, 'galaxy')}`;
    } else {
        // In the event that the page is accessed via an unexpected subdomain, use default backend
        console.log(`Unexpected subdomain, defaulting backend to ${galaxy_path}`);
    }
} else {
    // Dev only
    //galaxy_path = 'foobar';
    //galaxy_path = '/galaxy';
    //galaxy_path = `${window.location.protocol}//galaxy.dev.islandcompare.ca`;
    galaxy_path = `${window.location.protocol}//localhost:8000`;
}
