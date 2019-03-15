<!--template>

</template-->

<script>
    import HistoryItemFunction from '../HistoryItemFunction';
    export default {
        extends: HistoryItemFunction,
        name: "Remove",
        props: {
            icon: {
                type: String,
                default: 'icon-trash',
            },
            label: {
                type: String,
                default: 'Remove',
            },
            description: {
                type: String,
                default: 'Remove this item from the history',
            },
        },
    }
</script>

<style scoped>

</style>