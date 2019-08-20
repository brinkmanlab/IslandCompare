<!--template>

</template-->

<script>
    import FunctionIcon from "@/galaxy/misc/FunctionIcon";

    export default {
        extends: FunctionIcon,
        name: "ErrorInfo",
        props: {
            icon: {
                type: String,
                default: 'icon-errorinfo',
            },
            label: {
                type: String,
                default: 'Error Info',
            },
            description: {
                type: String,
                default: 'View error log',
            },
            action: {
                type: Function,
                default() {return ()=>{
                    this.item.model.get_error_log().then(log=>
                        this.$bvModal.msgBoxOk(this.$createElement('pre', {class:"error-log"}, [log]), {
                            title: 'Error',
                            size: 'xl',
                            buttonSize: 'sm',
                            okVariant: 'danger',
                            headerClass: 'p-2 border-bottom-0',
                            footerClass: 'p-2 border-top-0',
                            centered: true
                        })
                    );
                }}
            },
        },
    }
</script>

<style scoped>

</style>