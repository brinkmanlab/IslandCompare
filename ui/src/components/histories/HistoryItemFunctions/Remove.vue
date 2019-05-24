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
        methods: {
            onClick() {
                this.item.model.deleted = true;
                this.item.model.delete({url: this.item.model.history.contents_url});
            }
        },
    }
</script>

<style scoped>

</style>