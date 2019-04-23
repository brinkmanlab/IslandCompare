module.exports = {
    publicPath: '/islandcompare/',
    devServer: {
        disableHostCheck: true,
        proxy: {
            '^/api': {
                target: 'http://galaxy.brinkman.mbb.sfu.ca',
                changeOrigin: true
            },
            '^/user': {
                target: 'http://galaxy.brinkman.mbb.sfu.ca',
                changeOrigin: true
            }
        }
    },
}