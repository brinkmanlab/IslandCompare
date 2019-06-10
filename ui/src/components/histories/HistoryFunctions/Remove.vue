<!--template>

</template-->

<script>
    import HistoryFunction from '../HistoryFunction';
    export default {
        extends: HistoryFunction,
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
        },
        methods: {
            onClick() {
                this.model.delete({query: {purge: true}}); //TODO Remove purge=True?
                this.$emit('galaxy-history-deleted', this.model);
            }
        },
    }
</script>

<style scoped>

</style>