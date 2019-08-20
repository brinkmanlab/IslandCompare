<!--template>

</template-->

<script>
    import FunctionIcon from '../../misc/FunctionIcon';
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
            action: {
                type: Function,
                default() {
                    this.item.model.delete({query: {purge: true}}); //TODO Remove purge=True?
                    this.$emit('galaxy-history-deleted', this.item.model);
                },
            },
        },
    }
</script>

<style scoped>

</style>