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
        }),
        computed: {
        },
        methods: {
            update() {
                galaxy.workflows.WorkflowInvocation.$get({
                    params: {
                        url: this.model.workflow.url,
                        id: this.model.id,
                    },
                });
            }
        },
        mounted() {
            this.update();
        }
    }
</script>

<style scoped>

</style>