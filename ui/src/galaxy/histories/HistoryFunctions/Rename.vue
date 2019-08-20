<!--template>

</template-->

<script>
    import FunctionIcon from '../../misc/FunctionIcon';
    import History from "@/galaxy/histories/History";

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
                type: History,
                required: true,
            },
            action: {
                type: Function,
                default() {
                    this.$emit('galaxy-history-rename', this.item.model);
                },
            },
        },
    }
</script>

<style scoped>

</style>