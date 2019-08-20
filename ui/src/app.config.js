/*
Configuration options specific to the app
 */

export const workflow_name = "IslandCompare";
// TODO move to workflow input annotation: export const permitted_file_extensions = ['gbk', 'genbank', 'embl', 'gbff'];
export const base_path = '/islandcompare/'; // TODO change to subdomains so that base_path is root (/)
export const galaxy_path = '/galaxy'; // TODO enable CORS for galaxy and change this to domain name

//TODO move all global parameters here including workflow invocation data