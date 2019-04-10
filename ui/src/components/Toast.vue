<template>
    <div class="Toast" v-if="message" @click="hide"><span>{{ message }}</span></div>
</template>

<script>
    export default {
        name: "Toast",
        data: ()=>({
            timeoutHandle: null,
            message: "",
        }),
        methods: {
            show(message, timeout = 10000) {
                this.message = message;
                this.timeoutHandle = setTimeout(()=>this.message = "", timeout);
            },
            hide() {
                clearTimeout(this.timeoutHandle);
                this.message = "";
            }
        },
    }
</script>

<style scoped>
    .Toast {
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        border: red 2px dashed;
        background-color: lightcoral;
        padding-left: 3em;
        padding-right: 3em;
        width: fit-content;
        height: 4em;
    }
</style>