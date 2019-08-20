<!--template>

</template-->

<script>
    import FunctionIcon from '../../misc/FunctionIcon';
    import HistoryItem from "@/galaxy/histories/HistoryItem";
    export default {
        extends: FunctionIcon,
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
            item: {
                type: HistoryItem,
                required: true,
            },
            action: {
                type: Function,
                required: false,
                default() {
                    this.item.model.deleted = true;
                    this.item.model.delete({url: this.item.model.history.contents_url});
                    this.$emit('galaxy-history-item-deleted', this.item.model);
                },
            },
        },
    }
</script>

<style scoped>

</style>