<!--template>

</template-->

<script>
    import HistoryItemFunction from '../HistoryItemFunction';
    export default {
        extends: HistoryItemFunction,
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
        },
        methods: {
            onClick() {
                this.$emit('galaxy-history-item-rename', this.item.model);
            }
        },
    }
</script>

<style scoped>

</style>