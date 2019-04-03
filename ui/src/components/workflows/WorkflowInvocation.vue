<template>
    <div class="WorkflowInvocation">
        <span>{{ model.history.name }}</span>
        <span>{{ model.history.state }}</span>
        <time v-bind:datetime="model.history.update_time">{{ (new Date(model.history.update_time)).toLocaleDateString(undefined, { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric', hour: 'numeric', minute: 'numeric', second: 'numeric' }) }}</time>
        <slot v-bind:model="model"></slot>
        <div class="functions">
            <slot name="functions" v-bind="this" />
            <button @click="remove">Remove</button>
        </div>
    </div>
</template>

<script>
    import * as galaxy from "@/galaxy";

    export default {
        name: "WorkflowInvocation",
        props: {
            model: {
                type: galaxy.workflows.WorkflowInvocation,
                required: true,
            }
        },
        data() {return {
            pollHandle: null,
            stopPoll: ['ok','error'],
        }},
        computed: {
            //done() {
            //    return
            //}
        },
        methods: {
            remove() {
                galaxy.histories.History.$delete({params: {id: this.model.history.id}}).then(()=>{ //TODO Add purge=True?

                });
            }
        },
        mounted() {
            if (!this.model.steps.length) galaxy.workflows.WorkflowInvocation.$get({params:{url: this.model.workflow.url, id: this.model.id}});
            if (!this.stopPoll.includes(this.model.history.state)) {
                this.pollHandle = setInterval(()=>{
                    galaxy.histories.History.$get({params:{id: this.model.id}}).then(()=>{
                        if (this.stopPoll.includes(this.model.history.state)) {
                            clearInterval(this.pollHandle);
                            this.pollHandle = null;
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
    .functions {
        text-align: right;
    }
</style>