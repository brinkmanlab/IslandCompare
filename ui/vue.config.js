module.exports = {
    publicPath: '/islandcompare/',
    devServer: {
        //Documentation: https://github.com/chimurai/http-proxy-middleware#proxycontext-config
        disableHostCheck: true,
        proxy: {
            '^/galaxy/.*': {
                pathRewrite: (path, req)=>path.replace('/galaxy', ''), //eslint-disable-line
                target: 'http://galaxy.brinkman.mbb.sfu.ca/',
                changeOrigin: true,
                ws: false,
            },
        }
    },
}