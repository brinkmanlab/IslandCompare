<!--template>

</template-->

<script>
    import FunctionIcon from '../../misc/FunctionIcon';
    import History from "@/galaxy/histories/History";
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
                default: 'Remove this history',
            },
            item: {
                type: History,
                required: true,
            },
            action: {
                type: Function,
                default() {
                    this.item.model.deleted = true;
                    this.item.model.delete();
                    //this.$emit('galaxy-history-deleted', this.item.model); TODO 'this' is not the button Vue object
                },
            },
        },
    }
</script>

<style scoped>

</style>