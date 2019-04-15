<template>
    <ol class="WorkflowInvocationSteps">
        <li v-bind:is="WorkflowInvocationStep"
            v-for="step in models"
            v-bind:key="step.id"
            v-bind:model="step"
        ></li>
    </ol>
</template>

<script>
    //import * as galaxy from "@/galaxy";
    import WorkflowInvocationStep from './WorkflowInvocationStep'

    export default {
        name: "WorkflowInvocationSteps",
        components: [
            WorkflowInvocationStep
        ],
        props: {
            models: {
                type: Array,
                required: true,
            }
        }
    }
</script>

<style scoped>

</style>