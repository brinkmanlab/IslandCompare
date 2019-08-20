<!--template>

</template-->

<script>
    import FunctionIcon from "@/galaxy/misc/FunctionIcon";
    import HistoryItem from "@/galaxy/histories/HistoryItem";

    export default {
        extends: FunctionIcon,
        name: "Rename",
        props: {
            icon: {
                type: String,
                default: 'icon-rename',
            },
            label: {
                type: String,
                default: 'Rename',
            },
            description: {
                type: String,
                default: 'Rename this item',
            },
            item: {
                type: HistoryItem,
                required: true,
            },
            action: {
                type: Function,
                required: false,
                default() {
                    this.$emit('galaxy-history-item-rename', this.item.model);
                },
            },
        },
    }
</script>

<style scoped>

</style>