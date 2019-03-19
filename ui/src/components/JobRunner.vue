<template>
    <div class="JobRunner">
        <HistoryContents v-bind:history_id="current_history" filter="dataset" @operation-fail="e=>this.$emit('toast', e)"/>
        <input type="submit" />
    </div>
</template>

<script>
    import HistoryContents from './HistoryContents';
    export default {
        name: "JobRunner",
        components: {
            HistoryContents
        },
        data: ()=>{return {

        }},
        computed: {
            current_history() {
                return this.$store.state.galaxy.current_history
            }
        },
        methods: {

        }
    }
</script>

<style scoped>
    .JobRunner input[type=submit] {

    }
</style>