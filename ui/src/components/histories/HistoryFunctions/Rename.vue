<!--template>

</template-->

<script>
    import HistoryFunction from '../HistoryFunction';
    export default {
        extends: HistoryFunction,
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
                this.$emit('galaxy-history-rename', this.item.model);
            }
        },
    }
</script>

<style scoped>

</style>