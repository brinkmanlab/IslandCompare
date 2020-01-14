/*
Configuration options specific to the app
 */

import {publicPath} from "../vue.config";

export const workflow_name = "IslandCompare unpacked";
export const base_path = publicPath;
export let galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.ca`; // Default backend

const subdomain = /^islandcompare/;

if (process.env.NODE_ENV === 'production') {
    // Production only
    if (subdomain.test(window.location.hostname)) {
        galaxy_path = `${window.location.protocol}//${window.location.hostname.replace(subdomain, 'galaxy')}`;
    } else {
        // In the event that the page is accessed via an unexpected subdomain, use default backend
        console.log(`Unexpected subdomain, defaulting backend to ${galaxy_path}`);
    }
} else {
    // Dev only
    //galaxy_path = 'foobar';
    //galaxy_path = '/galaxy';
    galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.sfu.ca`;
    //galaxy_path = `${window.location.protocol}//gateway.cedar.computecanada.ca:8459`;
    //galaxy_path = `${window.location.protocol}//ec2-13-57-236-209.us-west-1.compute.amazonaws.com:81`;
}
