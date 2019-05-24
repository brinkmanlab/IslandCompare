<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="WorkflowInvocation">
        <History
                v-bind:model="model.history"
        >
            <template v-slot:functions="slot">
                <slot name="functions" v-bind="self"/>
            </template>
        </History>
    </div>
</template>

<script>
    import * as galaxy from "@/galaxy";
    import History from "../histories/History";
    export default {
        name: "WorkflowInvocation",
        components: {
            History,
        },
        props: {
            model: {
                type: galaxy.workflows.WorkflowInvocation,
                required: true,
            }
        },
        data() {return {
            self: this,
            end_states: ["error", "done"]
        }},
        methods: {
        },
        computed: {
            states() {
                return this.model.steps.reduce((acc, cur)=>{acc[cur.state] = (acc[cur.state] || 0) + 1; return acc}, {});
            },
            state() {
                //Touch properties for reaction system to detect
                //this.states.error;
                //this.states.new;
                //this.states.running;

                if (this.states.error) return "error";
                if (this.states.new) return "running";
                return "done";
            },
            done() {
                return this.state === "done" && this.outputs && Object.entries(this.outputs).length && Object.values(this.outputs).every(o => o.state === 'ok');
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
                            if (!hda.end_states.includes(hda.state)) {
                                hda.start_polling(()=>{
                                    return hda.end_states.includes(hda.state);
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
            if (!this.end_states.includes(this.state)) {
                this.model.start_polling(()=>{
                    if (this.end_states.includes(this.state)) {
                        this.$emit('workflow-completed', this);
                        return true;
                    }
                    return false;
                });
            }
        },
        beforeDestroy() {
            this.model.stop_polling();
        }
    }
</script>

<style scoped>

</style>