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
                default: 'icon-i-cursor',
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

            }
        },
    }
</script>

<style scoped>

</style>