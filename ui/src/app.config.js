/*
Configuration options specific to the app
 */

export const workflow_name = "IslandCompare unpacked";
export const base_path = '/';
//export const galaxy_path = '/galaxy';
//export const galaxy_path = `${window.location.protocol}//${window.location.hostname}:81`;
export const galaxy_path = `${window.location.protocol}//galaxy.pathogenomics.sfu.ca`;

if (process.env.NODE_ENV === 'production') {
    // Production only

} else {
    // Dev only

}