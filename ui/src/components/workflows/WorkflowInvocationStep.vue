<template>
    <li class="WorkflowInvocationStep">

    </li>
</template>

<script>
    import * as galaxy from "@/galaxy";
    export default {
        name: "WorkflowInvocationStep",
        props: {
            model: {
                type: galaxy.workflows.WorkflowInvocationStep,
                required: true,
            }
        }
    }
</script>

<style scoped>

</style>