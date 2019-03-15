<template>
    <div class="JobRunner">
        <HistoryList history_id="$galaxy.current_history"/>
        <input type="submit" />
    </div>
</template>

<script>
    import HistoryList from './HistoryList';
    export default {
        name: "JobRunner",
        components: {
            HistoryList
        },
        data: ()=>{return {

        }},
        methods: {

        }
    }
</script>

<style scoped>
    .JobRunner input[type=submit] {

    }
</style>