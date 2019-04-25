module.exports = {
    publicPath: '/islandcompare/',
    devServer: {
        //Documentation: https://github.com/chimurai/http-proxy-middleware#proxycontext-config
        disableHostCheck: true,
        proxy: {
            '^/(?!islandcompare).*': {
                target: 'http://galaxy.brinkman.mbb.sfu.ca/',
                changeOrigin: true,
            },
        }
    },
}