<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="galaxy-workflow-invocation">
        <History
                v-if="model.history"
                v-bind:model="model.history"
                @galaxy-history-deleted="model.stop_polling()"
        >
            <template v-slot="">
                <span class="galaxy-workflow-invocation-state">{{ state }}</span>
                <span class="galaxy-workflow-invocation-progress">
                    <b-progress class="w-100" v-bind:max="step_count()">
                        <b-progress-bar variant="success" v-bind:animated="false" v-bind:value="(states['scheduled'] || 0) + (states['ok'] || 0)">{{progress_label('ok')}}</b-progress-bar>
                        <b-progress-bar variant="success" v-bind:animated="true" v-bind:value="states['running']">{{progress_label('running')}}</b-progress-bar>
                        <b-progress-bar variant="danger" v-bind:animated="true" v-bind:value="states['error']">{{progress_label('error')}}</b-progress-bar>
                        <b-progress-bar variant="info" v-bind:animated="false" v-bind:value="(states['new'] || 0) + (states['queued'] || 0)">{{progress_label('new')}}</b-progress-bar>
                    </b-progress>
                </span>
            </template>
            <template v-slot:functions="slot">
                <slot name="functions" v-bind="self" />
                <ErrorInfo v-if="state === 'error'" v-bind:item="self" v-on="$listeners"/>
            </template>
        </History>
    </div>
</template>

<script>
    import * as galaxy from "@/galaxy/";
    import History from "../histories/History";
    import ErrorInfo from "@/galaxy/workflows/WorkflowInvocationFunctions/ErrorInfo";
    export default {
        name: "WorkflowInvocation",
        components: {
            History,
            ErrorInfo,
        },
        props: {
            model: {
                type: galaxy.workflows.WorkflowInvocation,
                required: true,
                //TODO validator to check that history and steps loaded
            },

        },
        data() {return {
            self: this,
        }},
        methods: {
            step_count() {
                return Object.values(this.states).reduce((a,b)=>a+b, 0);
            },
            progress_label(state) {
                if (Object.values(this.states).length === 0) return '';
                if (state === 'ok') {
                    if (this.done) return 'done';
                    if (state in this.states)
                        return this.states[state] + ' done';
                    else if ('scheduled' in this.states) return this.states['scheduled'] + ' done';
                }
                if (this.done) return '';
                if (state === 'running') return this.states[state] + ' running';
                if (state === 'new' && state in this.states) return this.states[state] + ' pending';
                if (state === 'new' && 'queued' in this.states) return this.states['queued'] + ' pending';
                if (state === 'error') return this.states[state] + ' failed';
            },
        },
        computed: {
            states() {
                return this.model.states();
            },
            state() {
                return this.model.aggregate_state();
            },
            done() {
                return this.outputs && Object.entries(this.outputs).length && Object.values(this.outputs).every(o => o.state === 'ok');
            },
        },
        asyncComputed: {
            outputs: {
                async get() {
                    let result = {};
                    if (this.model === null) return {};
                    for (let key of Object.keys(this.model.outputs)) {
                        let hda = galaxy.history_contents.HistoryDatasetAssociation.find(this.model.outputs[key].id);
                        if (hda) result[key] = hda;
                        else {
                            await galaxy.history_contents.HistoryDatasetAssociation.$get({params:{url: this.model.history.contents_url, id: this.model.outputs[key].id}});
                            hda = galaxy.history_contents.HistoryDatasetAssociation.find(this.model.outputs[key].id);
                            result[key] = hda;
                            if (!hda.constructor.end_states.includes(hda.state)) {
                                hda.start_polling(()=>{
                                    return hda.constructor.end_states.includes(hda.state);
                                });
                            }
                        }
                    }
                    return result;
                },
                default: {},
            },
        },
        mounted() {
            if (!this.model.constructor.end_states.includes(this.model.aggregate_state())) {
                this.model.start_polling(()=>{
                    if (this.model.constructor.end_states.includes(this.model.aggregate_state())) {
                        this.$emit('workflow-completed', this);
                        return true;
                    }
                    return false;
                }, {query: {view: "element", step_details: true}});
            }
        },
        beforeDestroy() {
            if (this.model)
                this.model.stop_polling();
        }
    }
</script>

<style scoped>
    .galaxy-workflow-invocation-progress {
        flex-wrap: nowrap;
    }
</style>