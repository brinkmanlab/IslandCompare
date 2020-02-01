<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <b-container class="Analysis" fluid>
        <b-row align-h="between">
            <b-col xl="6">
                <b-tabs v-model="current_tab" class="analysis-tabs">
                    <b-tab title="Recent Jobs">
                        <!-- TODO https://bootstrap-vue.js.org/docs/components/tabs#add-custom-content-to-tab-title -->
                        <!-- Shows running and completed jobs -->
                        <Jobs v-bind:invocations="invocations" ref="jobs" @click.native.prevent="$router.push('/history')">
                        </Jobs>
                    </b-tab>
                    <b-tab class="help" title="Instructions">
                        <p>
                            <a @click.prevent="start_tour('tour')" href="#" class="tutorial-start button-icon inline"><i class="icon icon-tutorial"></i> Tutorial</a>
                            Please see <b-link to="/about" >About</b-link> and <b-link to="/faq">FAQ</b-link> for more information.
                        </p>
                        <p>
                            Check out these <b-link :to="`visualize?src=${origin}/demo/listeria_sample_analysis.gff3`">example Listeria</b-link> or <b-link :to="`visualize?src=${origin}/demo/pseudomonas_sample_analysis.gff3`">example Pseudomonas</b-link> analyses.
                        </p>
                        <HTMLFragment :content="instructions" />
                        <p>
                            <a @click.prevent="show_api_key" href="#" class="show-api-key button-icon inline"><i class="icon icon-api"></i> Show API Key</a><i class="icon icon-info" title="The API key is used for the command line interface and other utilities that access the backend directly."></i>
                        </p>
                    </b-tab>
                </b-tabs>
            </b-col>
            <b-col xl="6">
                <JobRunner v-bind:history="history"
                           v-bind:workflow="workflow"
                           v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
                           @galaxy-workflow-invocation="invoked"
                />
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import { getConfiguredWorkflow, getUploadHistory, getInvocations, onInvocation, fetchStateAndUploadHistory } from "../app";
    // TODO async load these components
    import JobRunner from '../components/JobRunner';
    import Jobs from '../components/Jobs';

    import {updateRoute, api_key_key} from "../auth";
    import HTMLFragment from "../components/HTMLFragment";
    import instructions from "html-loader!@/assets/instructions.htm";

    export default {
        name: "Analysis",
        components: {
            HTMLFragment,
            JobRunner,
            Jobs,
        },
        data() { return {
            current_tab: 1,
            origin: window.location.origin,
            instructions,
            auth_fail: false,
        }},
        props: {
            tour: {
                type: String,
                default: '',
            }
        },
        methods: {
            start_tour(tour) {
                if (this.workflow && tour) this.$tours[tour].start();
            },
            invoked(invocation) {
                // Show jobs tab
                self.current_tab=0;
                // Add application tag to history
                onInvocation(invocation)
            },
            show_api_key() {
                // getAPIKey verifies key, get directly from store instead
                this.$bvModal.msgBoxOk(localStorage.getItem(api_key_key), {
                    title: 'API Key',
                    size: 'sm',
                    buttonSize: 'sm',
                    okVariant: 'info',
                    headerClass: 'p-2 border-bottom-0',
                    footerClass: 'p-2 border-top-0',
                    centered: true
                });
            },
            init(force = false) {
                if (this.auth_fail || force) {
                    this.auth_fail = false;
                    fetchStateAndUploadHistory(true).then(()=>{ // gidPromise can fail and then be reassigned by fetch. This happens late and is a race condition with below
                        updateRoute(this.$router, this.$route);
                    }).catch(() => { // TODO emit an event on $root and listen for it instead
                        this.auth_fail = true;
                    });
                }
            },
        },
        computed: {
            workflow: getConfiguredWorkflow,
            history: getUploadHistory,
            invocations() {
                const workflow = getConfiguredWorkflow();
                if (this.auth_fail) return [];
                if (!workflow || !workflow.invocationsFetched) return null;
                return getInvocations(workflow);
            },
        },
        watch: {
            invocations(last_val, new_val) {
                // Switch current tab to recent jobs if there are jobs
                if (new_val && new_val.length > 0 && (!last_val || last_val.length !== new_val.length)) {
                    this.current_tab = 0;
                }
            },
            workflow() {
                this.start_tour(this.tour);
            },
            tour() {
                if (this.tour) {
                    updateRoute(this.$router, this.$route);
                    if (this.workflow) this.$tours[this.tour].start();
                }
            }
        },
        activated() {
            this.init();
        },
        created() {
            this.init(true);
        },
        deactivated() {
            if (this.tour) this.$tours[this.tour].stop();
        }
    }
</script>

<style scoped>
    .analysis-tabs >>> .nav-tabs {
        font-size: 0.8em;
    }

    .analysis-tabs >>> .tab-content {
        border: 1px solid #dee2e6;
        border-top: 0;
        border-bottom-left-radius: 0.25rem;
        border-bottom-right-radius: 0.25rem;
        min-height: 70vh;
    }

    .JobRunner {
        padding-left: unset;
        padding-right: unset;
    }

    .JobRunner >>> .galaxy-workflow-parameter-dataset .b-table-sticky-header {
        min-height: 4rem;
    }

    .JobRunner >>> .galaxy-workflow-parameters .Input-datasets .b-table-sticky-header {
        height: 10rem;
    }

    .JobRunner >>> .input-group {
        font-size: 1em;
    }

    .JobRunner >>> .input-group * {
        font-size: inherit;
    }

    .help {
        /*max-width: 30em;*/
        padding: 1em;
        font-size: 0.9em;
    }

    .help >>> dl {
        counter-reset: instructions-counter;
    }

    .help >>> dt:before {
        content: counter(instructions-counter) '. ';
        counter-increment: instructions-counter;
    }

    .help .icon-info {
        font-size: 20px;
        color: var(--info);
    }

    .Jobs {
        padding: 0.5vw;
        border-collapse: separate;
    }

    .Jobs >>> thead {
        display: none;
    }

    .Jobs >>> .galaxy-workflow-invocation > :not(.galaxy-history) {
        display: none;
        /*visibility: collapse;*/
    }

    .Jobs >>> .galaxy-workflow-invocation > .galaxy-history > * {
        display: none;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-history-label, .Jobs >>> .galaxy-workflow-invocation .galaxy-workflow-invocation-progress:not(.new), .Jobs >>> .galaxy-workflow-invocation-state.new {
        display: table-cell;
        font-size: 0.7em;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-workflow-invocation-progress {
        width: 100%;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-history-label {
        padding-left: 1em;
        padding-right: 1em;
        white-space: nowrap;
    }
</style>
