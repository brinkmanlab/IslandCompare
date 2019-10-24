/*
Configuration options specific to the app
 */

export const workflow_name = "IslandCompare unpacked";
export const base_path = '/';
//export const galaxy_path = '/galaxy';
export const galaxy_path = `${window.location.protocol}//${window.location.hostname.replace(/^islandcompare/, 'galaxy')}`;
//export const galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.sfu.ca`;
//export const galaxy_path = `${window.location.protocol}//gateway.cedar.computecanada.ca:8459`;

if (process.env.NODE_ENV === 'production') {
    // Production only

} else {
    // Dev only

}