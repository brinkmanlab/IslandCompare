<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <JobManager v-if="user_id"
                v-bind:permitted_file_extensions="permitted_file_extensions"
                v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
    >
        <template v-slot:workflow_params="params">
            <label>Minimum island size<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
            <label>Minimum homologous region<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
            <label>Minimum cluster size<b-form-input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
        </template>
        <template v-slot:invocation_functions="slot">
            <template v-if="slot.done">
                <b-link v-bind:href="`/plugins/visualizations/islandcompare/show?dataset_id=${slot.outputs['IslandCompare Result'].id}` | auth">Visualize</b-link>
                <b-link v-bind:href="`/datasets/${slot.outputs['IslandCompare Result'].id}/display?to_ext=gff3` | auth">Download</b-link>
            </template>
        </template>
    </JobManager>
    <div v-else>
    </div>
</template>

<script>
    import Vue from 'vue'
    import JobManager from "@/components/JobManager";
    import uuid from 'uuid/v1';
    import axios from 'axios';

    import { permitted_file_extensions } from "@/app.config";

    export default {
        name: "Analysis",
        components: {
            JobManager
        },
        data() {return{
            user_id: null,
            permitted_file_extensions: permitted_file_extensions,
        }},
        asyncComputed: {

        },
        methods: {
            async uuid(){ // TODO move this into a separate js file
                // register filter to append ?uuid= to urls
                Vue.filter('auth', value=>value + (value.includes('?') ? '&' : '?') + 'uuid=' + id);
                if (this.user_id !== null) return this.user_id;
                let id = (new URLSearchParams(location.search)).get('uuid');
                if (!id) {
                    let m = document.cookie.match(/galaxysession_user_uuid=([^;]+)/);
                    if (m) {
                        //If uuid cookie set, use that.
                        id = m[1];
                    } else {
                    //else generate uuid
                    // TODO move this logic to backend
                    id = uuid();
                    let tag = location.search.lastIndexOf('#');
                    let url = "";
                    if (tag >= 0) {
                        url = location.search.slice(0, tag) + (location.search.includes('?') ? '&' : '?') + id + location.search.slice(tag);
                    } else {
                        url = location.search + (location.search.includes('?') ? '&' : '?') + "uuid=" + id;
                    }
                    history.replaceState(history.state, "Analysis", url);
                    alert("Be sure to bookmark this page to return to your work. The URL is unique to you."); //TODO replace with a html popup
                    }
                }
                document.cookie = `galaxysession_user_uuid=${id};path=/;max-age=31536000`;

                //this is a bandaid to get a session key from the galaxy frontend rather than an api key, api keys are not available to remote auth users
                let response = await axios.get('/user', { //eslint-disable-line
                    params: {
                        uuid: id,
                    }
                });

                /*/ fetch api key
                let response = await axios.get('/api/users', { //eslint-disable-line
                    params: {
                        uuid: id,
                    }
                });

                response = await axios.get(`/api/users/${response.data[0].id}/api_key/inputs`);
                */


                //CommonModel.methodConf.http.params.uuid = id;
                this.user_id = id;
                return id;
            }
        },
        mounted() {
            this.uuid()
        }
    }
</script>

<style scoped>

</style>