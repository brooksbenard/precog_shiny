Sidebar <-dashboardSidebar(
            sidebarMenu(
              id = "tabs",
              menuItem("Home", tabName = "home", icon = icon("home")
                       # menuSubItem("iPRECOG", tabName = "iPRECOG_home")
                       ),
              # menuItem("iPRECOG - Home", tabName = "iPRECOG_home", icon = icon("house-user")),
              menuItem("Data", tabName = "data", icon = icon("database")),
              menuItem("Analysis", tabName = "metaz", icon = icon("chart-line")),
                      
              menuItem("Individual Analysis", tabName = "individual_analysis", icon = icon("chart-pie")),
              menuItem("iPRECOG", tabName = "iPRECOG", icon = icon("syringe"),
                       menuSubItem("Home", tabName = "iPRECOG_home", icon = icon("home")),
                       menuSubItem('Signature Matrix', tabName = 'signature', icon = icon('file-signature')),
                       menuSubItem('iMetaZ matrix', tabName = 'iMetaz', icon = icon('chart-line')),
                       menuSubItem('Immune fractions', tabName = 'immune-fraction', icon = icon('disease'))
                       ),
              menuItem("Download", tabName = "downloads", icon = icon("cloud-download-alt")),
              # menuItem("Select Ressources", tabName = "ressources", icon = icon("mouse-pointer")),
              menuItem("Contact", tabName = "contact", icon = icon("address-book"))
              # menuItem('login', tabName = 'login', icon = icon('home'))
            )
          )
