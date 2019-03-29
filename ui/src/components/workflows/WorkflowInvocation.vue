<template>
    <div class="WorkflowInvocation">
        <span>{{ model.history.name }}</span>
        <span>{{ model.state }}</span>
        <time v-bind:datetime="model.update_time">{{model.update_time}}</time>
        <slot v-bind:model="model"></slot>
        <div class="functions">
            <slot name="functions" />
            <button>View</button>
            <button>X</button>
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
        }
    }
</script>

<style scoped>
    .functions {
        text-align: right;
    }
</style>