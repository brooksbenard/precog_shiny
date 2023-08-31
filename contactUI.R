source('code/utils.R')

# UI and Shiny Inputs
contactItem <- tabItem(
  
  tabName = "contact",
  fluidPage(

    includeCSS('code/css/contact.css'),
    
    h1("Contact", id = "contact"),
    
    fluidRow(
      box(
        # includeCSS('code/css/contact.css'),
        htmlOutput('andrew', inline = T),
        span(icon(name = 'envelope',id = 'env'), 'agentles at stanford dot edu'),
        title = h2(a('Andrew Gentles', href = 'http://med.stanford.edu/profiles/andrew-gentles/')),
      ),
      box(
        htmlOutput('aaron', inline = T),
        span(icon(name = 'envelope', id = 'env'), 'amnewman at stanford dot edu'),
        title = h2(a('Aaron Newman', href = 'http://med.stanford.edu/profiles/Aaron_Newman/'))
      ),
      box(
        htmlOutput('chih', inline = T),
        span(icon(name = 'envelope', id = 'env'), 'cll at stanford dot edu'),
        title = h2(a('Chih Long Liu', href = 'http://med.stanford.edu/profiles/Chih_Liu/'))
      ),
      box(
        htmlOutput('sylvia', inline = T),
        span(icon(name = 'envelope', id = 'env'), 'sylvia dot plevritis at stanford dot edu'),
        title = h2(a('Sylvia K. Plevritis', href = 'http://med.stanford.edu/profiles/sylvia-plevritis/'))
      ),
      
      box(

        htmlOutput('ash', inline = T),
        span(icon(name = 'envelope', id = 'env'), 'arasha at stanford dot edu'),
        title = h2(a('Ash Alizadeh', href = 'http://med.stanford.edu/profiles/Arash_Alizadeh/'))
      ),
      
    ),
    
    br(), 
    
    # includeHTML('code/html/contact.html')
    
    )
    
  )

# Shiny Outputs
andrew.output <- renderUI({HTML('<a title = "Dr Gentles" href = http://med.stanford.edu/profiles/andrew-gentles/><img id = "contact_img", src="andrew.jpg", height="100px, width = 100px"
                            style="vertical-align:middle;margin:0px 0px"/></a>')})
aaron.output <- renderUI({HTML('<img id = "contact_img", src="aaron.jpg", height="100px, width = 100px"
                            style="vertical-align:middle;margin:0px 0px"/>')})
chih.output <- renderUI({HTML('<img id = "contact_img", src="chih.jpg", height="100px, width = 100px"
                            style="vertical-align:middle;margin:0px 0px"/>')})
sylvia.output <- renderUI({HTML('<img id = "contact_img", src="sylvia.jpg", height="100px, width = 100px"
                            style="vertical-align:middle;margin:0px 0px"/>')})
ash.output <- renderUI({HTML('<img id = "contact_img", src="ash.jpg", height="100px, width = 100px"
                            style="vertical-align:middle;margin:0px 0px"/>')})
