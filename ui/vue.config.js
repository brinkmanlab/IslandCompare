module.exports = {
    publicPath: '/',
    devServer: {
        //Documentation: https://github.com/chimurai/http-proxy-middleware#proxycontext-config
        disableHostCheck: true,
        proxy: {
            '^/galaxy/.*': {
                //pathRewrite: (path, req)=>path.replace('/galaxy', ''), //eslint-disable-line
                //target: 'http://galaxy.brinkman.mbb.sfu.ca/',
                //target: 'https://gateway.cedar.computecanada.ca:8457/',
                target: 'http://ec2-13-57-236-209.us-west-1.compute.amazonaws.com/',
                changeOrigin: true,
                ws: false,
            },
        }
    },
}