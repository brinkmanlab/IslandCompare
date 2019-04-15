<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="WorkflowInvocation">
        <History
                v-bind:model="model.history"
                @history-completed="update"
        >
            <template v-slot:functions="slot">
                <slot name="functions" v-bind="model" />
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
        data: ()=>({
            self: this,
            pollHandle: null,
        }),
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
        },
        mounted() {
            if (!this.stopPoll.includes(this.model.state)) {
                this.pollHandle = setInterval(()=>{
                    galaxy.workflows.WorkflowInvocation.$get({
                        params: {
                            url: this.model.workflow.url,
                            id: this.model.id,
                        },
                    }).then(()=>{
                        if (this.model.end_states.includes(this.model.state)) {
                            clearInterval(this.pollHandle);
                            this.pollHandle = null;
                            this.$emit('history-completed', this);
                        }
                    });
                }, 10000);
            }
        },
        beforeDestroy() {
            if (this.pollHandle !== null) clearInterval(this.pollHandle);
        }
    }
</script>

<style scoped>

</style>