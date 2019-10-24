/*
Configuration options specific to the app
 */

export const workflow_name = "IslandCompare unpacked";
export const base_path = '/';
export let galaxy_path;

if (process.env.NODE_ENV === 'production') {
    // Production only
    galaxy_path = `${window.location.protocol}//${window.location.hostname.replace(/^islandcompare/, 'galaxy')}`;
} else {
    // Dev only
    //export const galaxy_path = '/galaxy';
    galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.sfu.ca`;
    //export const galaxy_path = `${window.location.protocol}//gateway.cedar.computecanada.ca:8459`;
}