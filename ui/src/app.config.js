/*
Configuration options specific to the app
 */

export const workflow_name = "IslandCompare unpacked";
export const base_path = '/';
//export const galaxy_path = '/galaxy'; // TODO enable CORS for galaxy and change this to domain name
//export const galaxy_path = `${window.location.protocol}//${window.location.hostname}:81`;
export const galaxy_path = 'http://ec2-13-57-236-209.us-west-1.compute.amazonaws.com:81';

if (process.env.NODE_ENV === 'production') {
    // Production only

} else {
    // Dev only

}