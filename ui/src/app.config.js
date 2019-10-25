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
    //galaxy_path = '/galaxy';
    //galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.sfu.ca`;
    galaxy_path = `${window.location.protocol}//gateway.cedar.computecanada.ca:8459`;
    //galaxy_path = `${window.location.protocol}//ec2-13-57-236-209.us-west-1.compute.amazonaws.com:81`;
}